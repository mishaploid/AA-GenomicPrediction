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
### Included first two principal components to account for population structure (--covar)
### Run all phenotypes simultaneously (--mpheno -1)
### --dentist YES will pad missing phenotype values

rule gblup_h2:
    input:
        pheno = config["pheno_file"],
        kins = "models/gblup/kinships.grm.id"
    output:
        "models/gblup_h2/{trait}.gblup.reml"
    params:
        out = "models/gblup_h2/{trait}.gblup",
        trait = "{trait}",
        grm = "models/gblup/kinships"
    run:
        shell("{ldak} --reml {params.out} \
        --pheno {input.pheno} \
        --pheno-name {params.trait} \
        --grm {params.grm}")

rule gblup_h2_pcs:
    input:
        pheno = "data/processed/pheno_file_pcs",
        kins = "models/gblup/kinships.grm.id"
    output:
        "models/gblup_pc50/{trait}.gblup.reml"
    params:
        out = "models/gblup_pc50/{trait}.gblup",
        trait = "{trait}",
        grm = "models/gblup/kinships"
    run:
        shell("{ldak} --reml {params.out} \
        --pheno {input.pheno} \
        --pheno-name {params.trait} \
        --grm {params.grm}")
        # if {wildcards.trait} == ('glu', 'gly', 'val', 'BCAA', 'gly_t', 'val_t'):
        #     print("Including PC1 and PC2 for {trait}")
        #     shell("{ldak} --reml {params.out} \
        #     --pheno {input.pheno} \
        #     --pheno-name {params.trait} \
        #     --grm {params.grm} \
        #     --covar {input.pc2}")
        # else:
        #     if {wildcards.trait} == 'met':
        #         print("Including PC1 {trait}")
        #         shell("{ldak} --reml {params.out} \
        #         --pheno {input.pheno} \
        #         --pheno-name {params.trait} \
        #         --grm {params.grm} \
        #         --covar {input.pc1}")
        #     else:
        #         print("Including no PCs for {trait}")


### GBLUP for genomic prediction - REML model
### includes cross-validation
### Included first two principal components to account for population structure (--covar)
### Run all phenotypes simultaneously (--mpheno -1)
###     When there are multiple wildcards in output name, need to escape other
###     wildcards (i.e. not the trait id) with {{}}
### --dentist YES will pad missing phenotype values

rule gblup_reml:
    input:
        pheno = config["pheno_file"],
        kins = "models/gblup/kinships.grm.id",
        pc1 = "data/processed/pca.1",
        pc2 = "data/processed/pca.2"
    output:
        reml = expand("models/gblup/{{trait}}.cv{cv}.{index}.reml", cv = CV, index = INDEX)
    params:
        out = "models/gblup/{trait}.cv",
        trait = "{trait}",
        keep = "data/processed/cross_validation/cv",
        grm = "models/gblup/kinships"
    run:
        for i in CV:
            for j in INDEX:
                shell("{ldak} --reml {params.out}{i}.{j} \
                --pheno {input.pheno} \
                --pheno-name {params.trait} \
                --grm {params.grm} \
                --keep {params.keep}{i}.train{j}")
                # if {wildcards.trait} == ('glu', 'gly', 'val', 'BCAA', 'gly_t', 'val_t'):
                #     print("Including PC1 and PC2 for {trait}")
                #     shell("{ldak} --reml {params.out}{i}.{j} \
                #     --pheno {input.pheno} \
                #     --pheno-name {params.trait} \
                #     --grm {params.grm} \
                #     --keep {params.keep}{i}.train{j} \
                #     --covar {input.pc2}")
                # else:
                #     if {wildcards.trait} == 'met':
                #         print("Including PC1 for {trait}")
                #         shell("{ldak} --reml {params.out}{i}.{j} \
                #         --pheno {input.pheno} \
                #         --pheno-name {params.trait} \
                #         --grm {params.grm} \
                #         --keep {params.keep}{i}.train{j} \
                #         --covar {input.pc1}")
                #     else:
                #         print("Including no PCs for {trait}")


### Step 4d: Calculate BLUPs
### extract BLUPs and calculate scores
### use output to check correlations with true values (predictive ability)

rule gblup_blup:
    input:
        pheno = "data/processed/pheno_file",
        model = expand("models/gblup/{{trait}}.cv{cv}.{index}.reml", cv = CV, index = INDEX)
    output:
        blup = expand("models/gblup/{{trait}}.cv{cv}.{index}.blup", cv = CV, index = INDEX),
        score = expand("models/gblup/{{trait}}.cv{cv}.{index}.profile", cv = CV, index = INDEX)
    params:
        out = "models/gblup/{trait}.cv",
        grm = "models/gblup/kinships",
        reml = "models/gblup/{trait}.cv",
        bfile = config["bfile"],
        keep = "data/processed/cross_validation/cv",
        trait = "{trait}",
        pc2 = "data/processed/pca.2"
    run:
        for i in CV:
            for j in INDEX:
                shell("{ldak} --calc-blups {params.out}{i}.{j} \
                --grm {params.grm} \
                --remlfile {params.reml}{i}.{j}.reml \
                --bfile {params.bfile}")
                shell("{ldak} --calc-scores {params.out}{i}.{j} \
                --bfile {params.bfile} \
                --power 0 \
                --scorefile {params.out}{i}.{j}.blup \
                --keep {params.keep}{i}.test{j} \
                --pheno {input.pheno} \
                --pheno-name {params.trait}")
                # if {wildcards.trait} == ('glu', 'gly', 'val', 'BCAA', 'gly_t', 'val_t', 'ile_AspFam_aspCorr'):
                #     print("Including PC1 and PC2 for {trait}")
                #     shell("{ldak} --calc-blups {params.out}{i}.{j} \
                #     --grm {params.grm} \
                #     --remlfile {params.reml}{i}.{j}.reml \
                #     --bfile {params.bfile} \
                #     --covar {params.pc2}")
                # else:
                #     if {wildcards.trait} == 'met':
                #         print("Including PC1 for {trait}")
                #         shell("{ldak} --calc-blups {params.out}{i}.{j} \
                #         --grm {params.grm} \
                #         --remlfile {params.reml}{i}.{j}.reml \
                #         --bfile {params.bfile} \
                #         --covar {params.pc1}")
                #     else:

rule gblup_results:
    input:
        expand("models/gblup/{trait}.cv5.10.blup", trait = TRAIT)
    output:
        "reports/gblup.RData"
    run:
        shell("Rscript src/03_summarize_gblup.R")
