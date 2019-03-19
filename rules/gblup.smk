# Run GBLUP model in LDAK

### calculate kinships
# use power 0 based on recommendations for prediction
# --ignore-weights YES to assign an effect to each marker

rule gblup_kins:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam"
    output:
        grm = "models/gblup/kinships.grm.id"
    params:
        bfile = config["bfile"],
        weights = "data/processed/sections/weights.short",
        outdir = "models/gblup/kinships"
    run:
        shell("{ldak} --calc-kins-direct {params.outdir} \
        --bfile {params.bfile} \
        --ignore-weights YES \
        --power 0")

### GBLUP model to estimate genetic variances
### NOTE: Inclusion of principal components based on previous model selection results
### --dentist YES will pad missing SNP values

rule gblup_h2:
    input:
        pheno = config["pheno_file"],
        kins = "models/gblup/kinships.grm.id",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    output:
        "models/gblup_h2/gblup_{trait}.reml"
    params:
        out = "models/gblup_h2/gblup_{trait}",
        trait = "{trait}",
        grm = "models/gblup/kinships",
    run:
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.out} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --grm {params.grm} \
            --covar {input.pc1}")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --covar {input.pc2}")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm}")

### GBLUP for genomic prediction - REML model
### includes cross-validation

rule gblup_reml:
    input:
        pheno = config["pheno_file"],
        kins = "models/gblup/kinships.grm.id",
        keep = "data/processed/cross_validation/cv_{cv}.train{index}"
    output:
        reml = "models/gblup/cv_{cv}_{index}_{trait}.reml"
    params:
        out = "models/gblup/cv_{cv}_{index}_{trait}",
        trait = "{trait}",
        grm = "models/gblup/kinships",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    run:
        if {wildcards.trait} == 13:
            print("Including PC1 for {trait}")
            shell("{ldak} --reml {params.out} \
            --pheno {input.pheno} \
            --mpheno {params.trait} \
            --grm {params.grm} \
            --keep {input.keep} \
            --covar {params.pc1}")
        else:
            if {wildcards.trait} == (5, 6, 18, 26, 43, 54):
                print("Including PC1 and PC2 for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --keep {input.keep} \
                --covar {params.pc2}")
            else:
                print("Including no PCs for {trait}")
                shell("{ldak} --reml {params.out} \
                --pheno {input.pheno} \
                --mpheno {params.trait} \
                --grm {params.grm} \
                --keep {input.keep}")

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
        bfile = config["bfile"],
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

rule gblup_results:
    output:
        "reports/gblup.RData"
    run:
        shell("Rscript src/03_summarize_gblup.R")
