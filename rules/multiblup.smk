# Run multiBLUP model for pathways

# export snplists for pathways
rule multiblup_pathways:
    input:
        "data/processed/gene_list_tair10.txt"
    output:
        "data/processed/pathways/{pathway}/list1"
    params:
        pathway = "{pathway}",
        bincode = bincode
    run:
        shell("Rscript src/05_select_pathway_snps.R {params.pathway} \"{params.bincode}\"")

# calculate kinships (power = 0, ignore weights YES)

rule multiblup_kins:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam",
        snps = "data/processed/pathways/{pathway}/list1"
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
### When there are multiple wildcards in output name, need to escape other
###   wildcards (i.e. not the trait id) with {{}}
### --dentist YES will pad missing phenotype values

rule multiblup_h2:
    input:
        pheno = config["pheno_file"],
        mgrm = expand("data/processed/pathways/{pathway}/partition.list", pathway = PATHWAYS.keys())
    output:
        expand("models/multiblup_h2/{pathway}/{{trait}}.multiblup.reml", pathway = PATHWAYS.keys())
    params:
        out = "models/multiblup_h2/",
        trait = "{trait}",
        mgrm = "data/processed/pathways/"
    run:
        for p in PATHWAYS.keys():
            shell("{ldak} --reml {params.out}{p}/{params.trait}.multiblup \
            --pheno {input.pheno} \
            --pheno-name {params.trait} \
            --mgrm {params.mgrm}{p}/partition.list \
            --constrain NO")

################################################################################
# MultiBLUP - genomic prediction

### MultiBLUP REML - predictive ability

rule multiblup_reml:
    input:
        pheno = config["pheno_file"],
        mgrm = expand("data/processed/pathways/{pathway}/partition.list", pathway = PATHWAYS.keys())
    output:
        reml = expand("models/multiblup/{pathway}/{{trait}}.cv5.10.reml", pathway = PATHWAYS.keys())
    params:
        out = "models/multiblup/",
        keep = "data/processed/cross_validation/cv",
        trait = "{trait}",
        mgrm = "data/processed/pathways/"
    run:
        for p in PATHWAYS.keys():
            for i in CV:
                for j in INDEX:
                    shell("{ldak} --reml {params.out}{p}/{params.trait}.cv{i}.{j} \
                    --pheno {input.pheno} \
                    --pheno-name {params.trait} \
                    --mgrm {params.mgrm}{p}/partition.list \
                    --keep {params.keep}{i}.train{j} \
                    --constrain NO")

### MultiBLUP - calculate BLUPs
### extract BLUPs from REML and calculate scores
### use output to check correlations with true values (predictive ability)

rule multiblup_blup:
    input:
        pheno = config["pheno_file"],
        model = expand("models/multiblup/{pathway}/{{trait}}.cv5.10.reml", pathway = PATHWAYS.keys()),
        check = expand("models/multiblup/{pathway}/{{trait}}.cv5.10.reml", pathway = PATHWAYS.keys())
    output:
        blup = expand("models/multiblup/{pathway}/{{trait}}.cv5.10.blup", pathway = PATHWAYS.keys()),
        score = expand("models/multiblup/{pathway}/{{trait}}.cv5.10.profile", pathway = PATHWAYS.keys())
    params:
        out = "models/multiblup/",
        mgrm = "data/processed/pathways/",
        bfile = config["bfile"],
        keep = "data/processed/cross_validation/cv",
        trait = "{trait}"
    run:
        for p in PATHWAYS.keys():
            for i in CV:
                for j in INDEX:
                    shell("{ldak} --calc-blups {params.out}{p}/{params.trait}.cv{i}.{j} \
                    --mgrm {params.mgrm}{p}/partition.list \
                    --remlfile {params.out}{p}/{params.trait}.cv{i}.{j}.reml \
                    --bfile {params.bfile}")
                    shell("{ldak} --calc-scores {params.out}{p}/{params.trait}.cv{i}.{j} \
                    --bfile {params.bfile} \
                    --power 0 \
                    --scorefile {params.out}{p}/{params.trait}.cv{i}.{j}.blup \
                    --keep {params.keep}{i}.test{j} \
                    --pheno {input.pheno} \
                    --pheno-name {params.trait}")

rule multiblup_results:
    input:
        expand("models/multiblup/{pathway}/{trait}.cv5.10.blup", pathway = PATHWAYS.keys(), trait = TRAIT)
    output:
        "reports/multiblup.RData"
    run:
        shell("Rscript src/06_summarize_multiblup.R")
