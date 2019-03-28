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
        pheno = config["pheno_file"],
        mgrm = "data/processed/pathways/{pathway}/partition.list"
    output:
        "models/multiblup_h2/{pathway}/{trait}.multiblup.reml"
    params:
        out = "models/multiblup_h2/{pathway}/{trait}.multiblup",
        trait = "{trait}",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        shell("{ldak} --reml {params.out} \
        --pheno {input.pheno} \
        --pheno-name {params.trait} \
        --mgrm {input.mgrm}")

        # if {wildcards.trait} == ('glu', 'gly', 'val', 'BCAA', 'gly_t', 'val_t'):
        #     print("Including PC1 and PC2 for {trait}")
        #     shell("{ldak} --reml {params.out_prefix} \
        #     --pheno {input.pheno} \
        #     --pheno-name {params.trait} \
        #     --mgrm {input.mgrm} \
        #     --covar {params.pc2}")
        # else:
        #     if {wildcards.trait} == 'met':
        #         print("Including PC1 for {trait}")
        #         shell("{ldak} --reml {params.out_prefix} \
        #         --pheno {input.pheno} \
        #         --pheno-name {params.trait} \
        #         --mgrm {input.mgrm} \
        #         --covar {params.pc1}")
        #     else:
        #         print("Including no PCs for {trait}")


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
        reml = expand("models/multiblup/{{pathway}}/{trait}.cv5.10.reml", trait = TRAIT)
    params:
        out = "models/multiblup/{pathway}/",
        keep = "data/processed/cross_validation/cv",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        for t in TRAIT:
            for i in CV:
                for j in INDEX:
                    shell("{ldak} --reml {params.out}{t}.cv{i}.{j} \
                    --pheno {input.pheno} \
                    --pheno-name {t} \
                    --mgrm {input.mgrm} \
                    --keep {params.keep}{i}.train{j}")

                # if {wildcards.trait} == ('glu', 'gly', 'val', 'BCAA', 'gly_t', 'val_t'):
                #     print("Including PC1 and PC2 for {trait}")
                #     shell("{ldak} --reml {params.out}{i}.{j} \
                #     --pheno {input.pheno} \
                #     --pheno-name {params.trait} \
                #     --mgrm {input.mgrm} \
                #     --covar {params.pc2} \
                #     --keep {params.keep}{i}.train{j}")
                # else:
                #     if {wildcards.trait} == 'met':
                #         print("Including PC1 and PC2 for {trait}")
                #         shell("{ldak} --reml {params.out}{i}.{j} \
                #         --pheno {input.pheno} \
                #         --pheno-name {params.trait} \
                #         --mgrm {input.mgrm} \
                #         --covar {params.pc1} \
                #         --keep {params.keep}{i}.train{j}")
                #     else:
                #         print("Including no PCs for {trait}")


### MultiBLUP - calculate BLUPs
### extract BLUPs from REML model and calculate scores
### use output to check correlations with true values (predictive ability)

rule multiblup_blup:
    input:
        pheno = "data/processed/pheno_file",
        model = expand("models/multiblup/{{pathway}}/{trait}.cv5.10.reml", trait = TRAIT),
        check = expand("models/multiblup/{{pathway}}/{trait}.cv5.10.reml", trait = TRAIT)
    output:
        blup = expand("models/multiblup/{{pathway}}/{trait}.cv5.10.blup", trait = TRAIT),
        score = expand("models/multiblup/{{pathway}}/{trait}.cv5.10.profile", trait = TRAIT)
    params:
        out = "models/multiblup/{pathway}/",
        mgrm = "data/processed/pathways/{pathway}/partition.list",
        bfile = config["bfile"],
        keep = "data/processed/cross_validation/cv"
    run:
        for t in TRAIT:
            for i in CV:
                for j in INDEX:
                    shell("{ldak} --calc-blups {params.out}{t}.cv{i}.{j} \
                    --mgrm {params.mgrm} \
                    --remlfile {params.out}{t}.cv{i}.{j}.reml \
                    --bfile {params.bfile}")
                    shell("{ldak} --calc-scores {params.out}{t}.cv{i}.{j} \
                    --bfile {params.bfile} \
                    --power 0 \
                    --scorefile {params.out}{t}.cv{i}.{j}.blup \
                    --keep {params.keep}{i}.test{j} \
                    --pheno {input.pheno} \
                    --pheno-name {t}")

rule multiblup_results:
    input:
        expand("models/multiblup/{pathway}/{trait}.cv5.10.blup", pathway = PATHWAYS.keys(), trait = TRAIT)
    output:
        "reports/multiblup.RData"
    run:
        shell("Rscript src/04_summarize_multiblup.R")
