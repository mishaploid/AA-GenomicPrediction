# Create training and test sets for cross validation
# used 10-fold cross validation with a one-fold hold out

# source dougspeed.com/wp-content/uploads/tartu_practical.pdf
# ./ldak5.XXX --cut-folds cv --bfile human --keep filter.keep --num-folds 10
# --cut-folds name of cross-validation test/training subsets
# --bfile input file basename (plink format)
# --keep use a subset of individuals
# --num-folds number of folds for cross-validation

rule get_cv:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam"
    output:
        test = "data/processed/cross_validation/cv{cv}.test{index}",
        train = "data/processed/cross_validation/cv{cv}.train{index}"
    params:
        bfile = config["bfile"],
        folds = "data/processed/cross_validation/cv{cv}",
        numfolds = 10
    run:
        shell("{ldak} --cut-folds {params.folds} \
        --bfile {params.bfile} \
        --num-folds {params.numfolds}")
