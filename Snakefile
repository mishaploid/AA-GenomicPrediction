# MultiBLUP workflow for free amino acid traits in Arabidopsis thaliana
# Sarah Turner-Hissong
# Last updated: 6 March 2019

# Phenotype data: Angelovici et al 2013 Plant Cell
# Genotype data: Li et al. 2010; Horton et al. 2012
# Pipeline follows MultiBLUP documentation at dougspeed.com/multiblup/

# CAVEATS (in no specific order)
# 1. Full workflow requires a substantial amount of disk space and is computationally intensive
#    Recommendation is to run just one or a few of the random subsets to start with, then scale up
# 2. Requires phenotype/genotype data to be previously downloaded
# 3. Paths to input files are not generic and will need to be updated if applied to different data
# 4. SNP sets need to be extracted beforehand
# 5. Principal components are included for specific traits (manual edits required to remove/modify)

# to execute workflow, use submit.sh file or run the following:
# snakemake --jobs N --rerun-incomplete --latency-wait 60 --cluster-config submit.json
# --cluster "sbatch -p {cluster.p} -o {cluster.o} --cpus-per-task {cluster.cpus-per-task}" -p &>> snakemake_logs/yymmdd.multiblup.log"

################################################################################
# define some variables

# list pathways that will be tested
# export SNP lists in advance using export_snplists.R
PATHWAYS = glob_wildcards("data/processed/pathways/{pathway}/list1").pathway
print(PATHWAYS)

# list range of traits (there are 54 traits in this analysis)
# number refers to a given column in the PLINK formatted phenotype file
TRAIT = list(range(1,55))

# list range of random subsets to include - we used 5000
# RANDOM = list(range(1,5001))
RANDOM = 1

# how many sets of cross validation - we used 5
# CV = list(range(1,6))
CV = 1
# fold for each cross validation - we used 10 fold (total of 50 cross validations)
# INDEX = list(range(1,11))
INDEX = 1

# directory for ldak software
ldak = "../software/ldak5.linux"

################################################################################
# this is a pseudo-rule that collects the target files
# to skip a step, can hash out the corresponding output here
rule all:
    input:
        # snp_weightings
        "data/processed/sections/weights.short",
        # pca
        "data/processed/pca.vect",
        # gblup_kins_h2
        "models/gblup/kinships.grm.id",
        # gblup_h2
        expand("models/gblup_h2/gblup_{trait}.reml", trait = TRAIT),
        # gblup_kins
        "models/gblup_h2/kinships.grm.id",
        # get_cv
        expand("data/processed/cross_validation/cv_{cv}.test{index}", \
        cv = CV, index = INDEX),
        expand("data/processed/cross_validation/cv_{cv}.train{index}", \
        cv = CV, index = INDEX),
        # gblup_reml
        expand("models/gblup/cv_{cv}_{index}_{trait}.reml", \
        cv = CV, index = list(range(1,11)), trait = TRAIT),
        # gblup_blup
        expand("models/gblup/cv_{cv}_{index}_{trait}.profile", \
        cv = CV, index = INDEX, trait = TRAIT),
        # calc_kins_h2
        expand("data/processed/pathways/{pathway}/partition.list", pathway = PATHWAYS),
        # multiblup_h2
        expand("models/multiblup_h2/{pathway}/multiblup_h2_{trait}.reml", \
        pathway = PATHWAYS, trait = TRAIT),
        # calc_kins_gp
        expand("models/multiblup_cv/{pathway}/partition.list", pathway = PATHWAYS),
        # multiblup_reml
        expand("models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}.reml", \
        pathway = PATHWAYS, cv = CV, index = INDEX, trait = TRAIT),
        # multiblup_blup
        expand("models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}.profile", \
        pathway = PATHWAYS, cv = CV, index = INDEX, trait = TRAIT),
        # calc_kins_control
        expand("data/processed/random_sets/c_{random}/partition.list", random = RANDOM),
        # reml_h2_control
        expand("models/reml_null/c_{random}/reml_h2_{trait}.reml", \
        random = RANDOM, trait = TRAIT)

################################################################################
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
        "data/processed/atregmap_clean.bed",
        "data/processed/atregmap_clean.bim",
        "data/processed/atregmap_clean.fam",
        "data/processed/atregmap_clean.map"
    params:
        indir = "data/raw/call_method_75_TAIR9",
        outdir = "data/processed/atregmap_clean"
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
        "data/processed/atregmap_clean.map"
    output:
        "data/processed/gene_list_tair10.txt",
        "data/processed/snp_gene_ids_tair10.txt"
    shell:
        "Rscript src/02_extract_gene_ids.R {input}"

################################################################################
# Step 1: calculate SNP weightings
# source: dougspeed.com/get-weightings/
# input: plink .bed/.bim/.fam files for all genomic SNPs
# output: partitions/weights.all and partitions/weights.short
# this is for the union of all SNPs (aka all genomic SNPs)
# note: both union of subsets and weightings include all genomic SNPs

rule snp_weightings:
    input:
        bed = "data/processed/input_nomissing.bed",
        bim = "data/processed/input_nomissing.bim",
        fam = "data/processed/input_nomissing.fam"
    output:
        "data/processed/sections/weights.all",
        "data/processed/sections/weights.predictors",
        "data/processed/sections/weights.short"
    params:
        bfile = "data/processed/input_nomissing"
    run:
        shell("{ldak} --cut-weights sections \
        --bfile {params.bfile} \
        --section-kb 0.25")
        shell("{ldak} --calc-weights-all sections \
        --bfile {params.bfile}")

################################################################################
# Step 2: run principal component analysis (PCA)
# use as covariates in reml model to account for population structure
# NOTE:
# before running script: awk < {input.bim} '{print $2}' > data/processed/all.snps
# after completion: awk '{print $1, $2, $3}' ../data/processed/pca.vect > ../data/processed/pca.1

rule pca:
    input:
        bed = "data/processed/input_nomissing.bed",
        bim = "data/processed/input_nomissing.bim",
        fam = "data/processed/input_nomissing.fam"
    output:
        pca = "data/processed/pca.vect",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2",
        pc3 = "data/processed/pca.3"
    params:
        bfile = "data/processed/input_nomissing",
        pruned = "data/processed/prune",
        outdir = "data/processed/pca"
    run:
        shell("plink --bfile {params.bfile} \
        --extract data/processed/all.snps \
        --make-founders require-2-missing \
        --indep-pairwise 200 10 .05 \
        --out {params.pruned}")
        shell("{ldak} --calc-kins-direct {params.outdir} \
        --bfile {params.bfile} \
        --extract {params.pruned}.prune.in \
        --ignore-weights YES \
        --power -0.25")
        shell("{ldak} --pca {params.outdir} \
        --grm {params.outdir}")
        shell("{ldak} --calc-pca-loads {params.outdir} \
        --pcastem {params.outdir} \
        --grm {params.outdir} \
        --bfile {params.bfile}")
        shell("awk '{print $1, $2, $3}' {output.pca} > {output.pc1}")
        shell("awk '{print $1, $2, $4}' {output.pca} > {output.pc2}")
        shell("awk '{print $1, $2, $5}' {output.pca} > {output.pc3}")

################################################################################
# Step 4: GBLUP model
### Step 4a: Estimate kinships for prediction using GBLUP
### Ignores SNP weightings and power set to 0

rule gblup_kins:
    input:
        bed = "data/processed/input_nomissing.bed",
        bim = "data/processed/input_nomissing.bim",
        fam = "data/processed/input_nomissing.fam"
    output:
        grm = "models/gblup/kinships.grm.id"
    params:
        bfile = "data/processed/input_nomissing",
        weights = "data/processed/sections/weights.short",
        outdir = "models/gblup/kinships"
    run:
        shell("{ldak} --calc-kins-direct {params.outdir} \
        --bfile {params.bfile} \
        --ignore-weights YES \
        --power 0")

### Step 3b: GBLUP model to estimate genetic variances
### NOTE: Inclusion of principal components based on previous model selection results

rule gblup_h2:
    input:
        pheno = "data/processed/pheno_file",
        kins = "models/gblup/kinships.grm.id"
    output:
        h2 = "models/gblup_h2/gblup_{trait}.reml"
    params:
        out = "models/gblup_h2/gblup_{trait}",
        trait = "{trait}",
        grm = "models/gblup_h2/kinships",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.out} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --grm {params.grm} \
            --covar {params.pc1} \
            --dentist YES")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --covar {params.pc2} \
                --dentist YES")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --dentist YES")

################################################################################
# Step 4: Genomic prediction (GBLUP model)

### Step 4b: Create training and test sets for cross validation
### used 10-fold cross validation with a one-fold hold out

# source dougspeed.com/wp-content/uploads/tartu_practical.pdf
# ./ldak5.XXX --cut-folds cv --bfile human --keep filter.keep --num-folds 10
# --cut-folds name of cross-validation test/training subsets
# --bfile input file basename (plink format)
# --keep use a subset of individuals
# --num-folds number of folds for cross-validation

rule get_cv:
    input:
        bed = "data/processed/input_nomissing.bed",
        bim = "data/processed/input_nomissing.bim",
        fam = "data/processed/input_nomissing.fam"
    output:
        test = "data/processed/cross_validation/cv_{cv}.test{index}",
        train = "data/processed/cross_validation/cv_{cv}.train{index}"
    params:
        bfile = "data/processed/input_nomissing",
        folds = "data/processed/cross_validation/cv_{cv}",
        numfolds = 10
    run:
        shell("{ldak} --cut-folds {params.folds} \
        --bfile {params.bfile} \
        --num-folds {params.numfolds}")

### Step 4c: REML - predictive ability
### run reml model for predictive ability

rule gblup_reml:
    input:
        pheno = "data/processed/pheno_file",
        kins = "models/gblup/kinships.grm.id"
    output:
        reml = "models/gblup/cv_{cv}_{index}_{trait}.reml"
    params:
        out = "models/gblup/cv_{cv}_{index}_{trait}",
        trait = "{trait}",
        grm = "models/gblup/kinships",
        keep = "data/processed/cross_validation/cv_{cv}.train{index}",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.out} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --grm {params.grm} \
            --keep {params.keep} \
            --covar {params.pc1}")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --keep {params.keep} \
                --covar {params.pc2}")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --keep {params.keep}")

### Step 4d: Calculate BLUPs
### extract BLUPs and calculate scores
### use output to check correlations with true values (predictive ability)

rule gblup_blup:
    input:
        pheno = "data/processed/pheno_file",
        reml = "models/gblup/cv_{cv}_{index}_{trait}.reml"
    output:
        blup = "models/gblup/cv_{cv}_{index}_{trait}.blup",
        score = "models/gblup/cv_{cv}_{index}_{trait}.profile"
    params:
        out = "models/gblup/cv_{cv}_{index}_{trait}",
        grm = "models/gblup/kinships",
        bfile = "data/processed/input_nomissing",
        keep = "data/processed/cross_validation/cv_{cv}.test{index}",
        trait = "{trait}"
    run:
        shell("{ldak} --calc-blups {params.out} \
        --grm {params.grm} \
        --remlfile {input.reml} \
        --bfile {params.bfile}")
        shell("{ldak} --calc-scores {params.out} \
        --bfile {params.bfile} \
        --power 0 \
        --scorefile {params.out}.blup \
        --keep {params.keep} \
        --pheno {input.pheno} \
        --mpheno {params.trait}")

################################################################################
# Step 5: MultiBLUP model

rule multiblup_kins:
    input:
        bed = "data/processed/input_nomissing.bed",
        bim = "data/processed/input_nomissing.bim",
        fam = "data/processed/input_nomissing.fam"
    output:
        list = "models/multiblup/{pathway}/partition.list",
        k1 = "models/multiblup/{pathway}/kinships.1.grm.details",
        k2 = "models/multiblup/{pathway}/kinships.2.grm.details"
    params:
        bfile = "data/processed/input_nomissing",
        outdir = "models/multiblup/{pathway}",
        prefix = "data/processed/pathways/{pathway}/list"
    run:
        # partition kinship matrix
        # --ignore-weights YES & --power 0 based on LDAK recommendations for prediction
        shell("{ldak} --cut-kins {params.outdir} \
        --bfile {params.bfile} \
        --partition-number 2 \
        --partition-prefix {params.prefix}")
        shell("{ldak} --calc-kins {params.outdir} \
        --bfile {params.bfile} \
        --partition 1 \
        --ignore-weights YES \
        --power 0")
        shell("{ldak} --calc-kins {params.outdir} \
        --bfile {params.bfile} \
        --partition 2 \
        --ignore-weights YES \
        --power 0")

### Step 5b: MultiBLUP - REML model to estimate genetic variances (heritability)

rule multiblup_h2:
    input:
        pheno = "data/processed/pheno_file",
        mgrm = "models/multiblup/{pathway}/partition.list"
    output:
        "models/multiblup_h2/{pathway}/multiblup_h2_{trait}.reml"
    params:
        out_prefix = "models/multiblup_h2/{pathway}/multiblup_h2_{trait}",
        trait = "{trait}",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.out_prefix} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --mgrm {input.mgrm} \
            --covar {params.pc1}")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.out_prefix} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --mgrm {input.mgrm} \
                --covar {params.pc2}")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.out_prefix} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --mgrm {input.mgrm}")

################################################################################
# Step 6: MultiBLUP - genomic prediction

### Step 6a: MultiBLUP kinship matrices - genomic prediction
### calculate kinship matrices for prediction model using pathway SNPs
### list1 = markers in the feature set
### list2 = remaining genomic SNPs not in the feature set


### Step 6b: MultiBLUP REML - predictive ability
### run reml model for predictive ability

rule multiblup_reml:
    input:
        pheno = "data/processed/pheno_file",
        mgrm = "models/multiblup_cv/{pathway}/partition.list"
    output:
        reml = "models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}.reml"
    params:
        out_prefix = "models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}",
        trait = "{trait}",
        keep = "data/processed/cross_validation/cv_{cv}.train{index}"
    run:
        # shell("../software/ldak5.linux --reml {params.out_prefix} --pheno {input.pheno} --mpheno {params.trait} --mgrm {input.mgrm} --keep {params.keep}")
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.out_prefix} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --mgrm {input.mgrm} \
            --covar {params.pc1} \
            --keep {params.keep}")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.out_prefix} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --mgrm {input.mgrm} \
                --covar {params.pc2} \
                --keep {params.keep}")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.out_prefix} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --mgrm {input.mgrm} \
                --keep {params.keep}")

### Step 6c: MultiBLUP - calculate BLUPs
### extract BLUPs from REML model and calculate scores
### use output to check correlations with true values (predictive ability)

rule multiblup_blup:
    input:
        pheno = "data/processed/pheno_file",
        reml = "models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}.reml"
    output:
        blup = "models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}.blup",
        score = "models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}.profile"
    params:
        out = "models/multiblup_cv/{pathway}/cv_{cv}_{index}_{trait}",
        mgrm = "models/multiblup_cv/{pathway}/partition.list",
        bfile = "data/processed/input_nomissing",
        keep = "data/processed/cross_validation/cv_{cv}.test{index}",
        trait = "{trait}"
    run:
        shell("{ldak} --calc-blups {params.out} \
        --mgrm {params.mgrm} \
        --remlfile {input.reml} \
        --bfile {params.bfile}")
        shell("{ldak} --calc-scores {params.out} \
        --bfile {params.bfile} \
        --power 0 \
        --scorefile {output.blup} \
        --keep {params.keep} \
        --pheno {input.pheno} \
        --mpheno {params.trait}")

################################################################################
# Step 7: Create an empirical null distribution
# AKA the most time consuming step...
# Random SNP sets need to be generated in advance and stored in separate folders

### generate a uniform distribution of SNP sizes
rule null_sampling:
    output:
        "data/interim/null_group_sizes.txt"
    script:
        "src/03_null_group_sizes.R"

### Step 7a: calculate kinships for control sets
### Only interested in partitioning variance
### Use weights and power -0.25

rule calc_kins_control:
    input:
        bed = "data/processed/input_nomissing.bed",
        bim = "data/processed/input_nomissing.bim",
        fam = "data/processed/input_nomissing.fam"
    output:
        list = "data/processed/random_sets/c_{random}/partition.list",
        k1 = "data/processed/random_sets/c_{random}/kinships.1.grm.details",
        k2 = "data/processed/random_sets/c_{random}/kinships.2.grm.details"
    params:
        bfile = "data/processed/input_nomissing",
        outdir = "data/processed/random_sets/c_{random}",
        prefix = "data/processed/random_sets/c_{random}/list",
        weights = "data/processed/sections/weights.short"
    run:
        shell("{ldak} --cut-kins {params.outdir} \
        --bfile {params.bfile} \
        --partition-number 2 \
        --partition-prefix {params.prefix}")
        shell("{ldak} --calc-kins {params.outdir} \
        --bfile {params.bfile} \
        --partition 1 \
        --weights {params.weights} \
        --power 0")
        shell("{ldak} --calc-kins {params.outdir} \
        --bfile {params.bfile} \
        --partition 2 \
        --weights {params.weights} \
        --power 0")

### Step 7b: REML model for control sets

rule null_h2:
    input:
        pheno = "data/processed/pheno_file",
        mgrm = "data/processed/random_sets/c_{random}/partition.list",
    output:
        out = "models/null_h2/c_{random}/reml_h2_{trait}.reml",
    params:
        prefix = "models/null_h2/c_{random}/reml_h2_{trait}",
        trait = "{trait}",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.prefix} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --mgrm {input.mgrm} \
            --dentist YES \
            --covar {params.pc1}")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.prefix} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --mgrm {input.mgrm} \
                --dentist YES \
                --covar {params.pc2}")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.prefix} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --mgrm {input.mgrm} \
                --dentist YES")

################################################################################
# THE END!
# Artwork: Clover patch by Joan G. Stark
#         ,,,                      ,,,
#        {{{}}    ,,,             {{{}}    ,,,
#     ,,, ~Y~    {{{}},,,      ,,, ~Y~    {{{}},,,
#    {{}}} |/,,,  ~Y~{{}}}    {{}}} |/,,,  ~Y~{{}}}
#     ~Y~ \|{{}}}/\|/ ~Y~  ,,, ~Y~ \|{{}}}/\|/ ~Y~  ,,,
#     \|/ \|/~Y~  \|,,,|/ {{}}}\|/ \|/~Y~  \|,,,|/ {{}}}
#     \|/ \|/\|/  \{{{}}/  ~Y~ \|/ \|/\|/  \{{{}}/  ~Y~
#     \|/\\|/\|/ \\|~Y~//  \|/ \|/\\|/\|/ \\|~Y~//  \|/
#     \|//\|/\|/,\\|/|/|// \|/ \|//\|/\|/,\\|/|/|// \|/
# jgs^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
