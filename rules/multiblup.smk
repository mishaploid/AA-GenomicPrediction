# Run multiBLUP model for pathways

# export snplists for pathways
rule multiblup_pathways:
    input:
        "data/processed/gene_list_tair10.txt"
    output:
        "data/processed/pathways/{pathway}/list1"
    params:
        "{pathway}",
        bincode
    run:
        shell("Rscript src/04_select_pathway_snps.R {params}")

# calculate kinships (power = 0, ignore weights YES)

rule multiblup_kins:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam"
    output:
        list = "data/processed/pathways/{pathway}/partition.list",
        k1 = "data/processed/pathways/{pathway}/kinships.1.grm.details",
        k2 = "data/processed/pathways/{pathway}/kinships.2.grm.details"
    params:
        bfile = config["bfile"],
        outdir = "data/processed/pathways/{pathway}",
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

### MultiBLUP - REML model to estimate genetic variances (heritability)
### Included first two principal components to account for population structure (--covar)
### Run all phenotypes simultaneously (--mpheno -1)
###     When there are multiple wildcards in output name, need to escape other
###     wildcards (i.e. not the trait id) with {{}}
### --dentist YES will pad missing phenotype values

rule multiblup_h2:
    input:
        pheno = "data/processed/pheno_file",
        mgrm = "data/processed/pathways/{pathway}/partition.list"
    output:
        expand("models/multiblup_h2/{{pathway}}/multiblup_h2.{trait}.reml", trait = TRAIT)
    params:
        out_prefix = "models/multiblup_h2/{pathway}/multiblup_h2",
        trait = "-1",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        shell("{ldak} --reml {params.out_prefix} \
        --pheno {input.pheno} \
        --mpheno {params.trait} \
        --mgrm {input.mgrm} \
        --covar {params.pc2} \
        --dentist YES")
        # if {wildcards.trait} == 13:
        #     print("Including PC1 for {trait}")
        #     shell("{ldak} --reml {params.out_prefix} \
        #     --pheno {input.pheno} \
        #     --mpheno {params.trait} \
        #     --mgrm {input.mgrm} \
        #     --covar {params.pc1}")
        # else:
        #     if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
        #         print("Including PC1 and PC2 for {trait}")
        #         shell("{ldak} --reml {params.out_prefix} \
        #         --pheno {input.pheno} \
        #         --mpheno {params.trait} \
        #         --mgrm {input.mgrm} \
        #         --covar {params.pc2}")
        #     else:
        #         print("Including no PCs for {trait}")
        #         shell("{ldak} --reml {params.out_prefix} \
        #         --pheno {input.pheno} \
        #         --mpheno {params.trait} \
        #         --mgrm {input.mgrm}")

################################################################################
# MultiBLUP - genomic prediction

### MultiBLUP REML - predictive ability
### Included first two principal components to account for population structure (--covar)
### Run all phenotypes simultaneously (--mpheno -1)
###     When there are multiple wildcards in output name, need to escape other
###     wildcards (i.e. not the trait id) with {{}}
### --dentist YES will pad missing phenotype values

rule multiblup_reml:
    input:
        pheno = "data/processed/pheno_file",
        mgrm = "data/processed/pathways/{pathway}/partition.list"
    output:
        reml = expand("models/multiblup/{{pathway}}/cv_{{cv}}_{{index}}.{trait}.reml", trait = TRAIT)
    params:
        out_prefix = "models/multiblup/{pathway}/cv_{cv}_{index}",
        trait = "-1",
        keep = "data/processed/cross_validation/cv_{cv}.train{index}",
        pc2 = "data/processed/pca.2"
    run:
        shell("{ldak} --reml {params.out_prefix} \
        --pheno {input.pheno} \
        --mpheno {params.trait} \
        --mgrm {input.mgrm} \
        --covar {params.pc2} \
        --keep {params.keep} \
        --dentist YES")
        # if {wildcards.trait} == 13:
        #     print("Including PC1 for {trait}")
        #     shell("{ldak} --reml {params.out_prefix} \
        #     --pheno {input.pheno} \
        #     --mpheno {params.trait} \
        #     --mgrm {input.mgrm} \
        #     --covar {params.pc1} \
        #     --keep {params.keep}")
        # else:
        #     if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
        #         print("Including PC1 and PC2 for {trait}")
        #         shell("{ldak} --reml {params.out_prefix} \
        #         --pheno {input.pheno} \
        #         --mpheno {params.trait} \
        #         --mgrm {input.mgrm} \
        #         --covar {params.pc2} \
        #         --keep {params.keep}")
        #     else:
        #         print("Including no PCs for {trait}")
        #         shell("{ldak} --reml {params.out_prefix} \
        #         --pheno {input.pheno} \
        #         --mpheno {params.trait} \
        #         --mgrm {input.mgrm} \
        #         --keep {params.keep}")

### MultiBLUP - calculate BLUPs
### extract BLUPs from REML model and calculate scores
### use output to check correlations with true values (predictive ability)

rule multiblup_blup:
    input:
        pheno = "data/processed/pheno_file",
        reml = "models/multiblup/{pathway}/cv_{cv}_{index}.{trait}.reml"
    output:
        blup = "models/multiblup/{pathway}/cv_{cv}_{index}.{trait}.blup",
        score = "models/multiblup/{pathway}/cv_{cv}_{index}.{trait}.profile"
    params:
        out = "models/multiblup/{pathway}/cv_{cv}_{index}.{trait}",
        mgrm = "data/processed/pathways/{pathway}/partition.list",
        bfile = config["bfile"],
        keep = "data/processed/cross_validation/cv_{cv}.test{index}",
        trait = "{trait}",
        pc2 = "data/processed/pca.2"
    run:
        shell("{ldak} --calc-blups {params.out} \
        --mgrm {params.mgrm} \
        --remlfile {input.reml} \
        --bfile {params.bfile} \
        --covar {params.pc2}")
        shell("{ldak} --calc-scores {params.out} \
        --bfile {params.bfile} \
        --power 0 \
        --scorefile {output.blup} \
        --keep {params.keep} \
        --pheno {input.pheno} \
        --mpheno {params.trait}")