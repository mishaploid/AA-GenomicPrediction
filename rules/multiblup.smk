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
