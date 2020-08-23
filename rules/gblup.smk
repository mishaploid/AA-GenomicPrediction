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
        --power 0 \
        --kinship-raw YES")

### GBLUP model to estimate genetic variances
### Run all phenotypes simultaneously (--mpheno -1)
### --dentist YES will pad missing phenotype values

rule gblup_h2:
    input:
        pheno = config["pheno_file"],
        kins = "models/gblup/kinships.grm.id",
        covar = config["covar"]
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
        --grm {params.grm} \
        --covar {input.covar} \
        --constrain YES \
        --reml-iter 500")

### GBLUP for genomic prediction - REML
### includes cross-validation
### When there are multiple wildcards in output name, need to escape other
###   wildcards (i.e. not the trait id) with {{}}
### --dentist YES will pad missing phenotype values

rule gblup_reml:
    input:
        pheno = config["pheno_file"],
        kins = "models/gblup/kinships.grm.id",
        covar = config["covar"]
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
                --keep {params.keep}{i}.train{j} \
                --covar {input.covar} \
                --constrain YES \
                --reml-iter 500")

### Step 4d: Calculate BLUPs
### extract BLUPs and calculate scores
### use output to check correlations with true values (predictive ability)

rule gblup_blup:
    input:
        pheno = config["pheno_file"],
        model = expand("models/gblup/{{trait}}.cv{cv}.{index}.reml", cv = CV, index = INDEX),
        covar = config["covar"]
    output:
        blup = expand("models/gblup/{{trait}}.cv{cv}.{index}.blup", cv = CV, index = INDEX),
        score = expand("models/gblup/{{trait}}.cv{cv}.{index}.profile", cv = CV, index = INDEX)
    params:
        out = "models/gblup/{trait}.cv",
        grm = "models/gblup/kinships",
        reml = "models/gblup/{trait}.cv",
        bfile = config["bfile"],
        keep = "data/processed/cross_validation/cv",
        trait = "{trait}"
    run:
        for i in CV:
            for j in INDEX:
                shell("{ldak} --calc-blups {params.out}{i}.{j} \
                --grm {params.grm} \
                --remlfile {params.reml}{i}.{j}.reml \
                --bfile {params.bfile} \
                --covar {input.covar}")
                shell("{ldak} --calc-scores {params.out}{i}.{j} \
                --bfile {params.bfile} \
                --power 0 \
                --scorefile {params.out}{i}.{j}.blup \
                --keep {params.keep}{i}.test{j} \
                --pheno {input.pheno} \
                --pheno-name {params.trait}")

rule gblup_results:
    input:
        blups = expand("models/gblup/{trait}.cv5.10.blup", trait = TRAIT),
        blues = config["pheno_file"]
    output:
        "reports/gblup.RData"
    run:
        shell("Rscript src/04_summarize_gblup.R {input.blues}")
