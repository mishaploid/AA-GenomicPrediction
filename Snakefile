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

# to execute workflow, use submit.sh file or run the following:
# snakemake --jobs N --rerun-incomplete --latency-wait 60 --cluster-config submit.json
# --cluster "sbatch -p {cluster.p} -o {cluster.o} --cpus-per-task {cluster.cpus-per-task}" -p &>> snakemake_logs/yymmdd.multiblup.log"

################################################################################
import numpy as np

# define some variables

# dictionary of pathways and MapMan bincodes that will be tested
PATHWAYS = {
    "protein_aa_activation": "29.1",
    "protein_synthesis": "29.2",
    "protein_targeting": "29.3",
    "protein_postrans": "29.4",
#    "protein_degradation": "29.5",
    "degradation_subtilases": "29.5.1",
    "degradation_autophagy": "29.5.2",
    "cys_protease": "29.5.3",
    "asp_protease": "29.5.4",
    "ser_protease": "29.5.5",
    "metallo_protease": "29.5.7",
    "degradation_AAA": "29.5.9",
    "degradation_ubiquitin": "29.5.11",
    "protein_folding": "29.6",
    "protein_glyco": "29.7",
    "protein_assembly": "29.8",
    "aa_transport": "34.3",
    "tca_cycle": "8",
    "e_alt_oxidase": "9.4",
    "isoprenoids": "16.1",
    "phenylpropanoids": "16.2",
    "N_containing": "16.4",
    "S_containing": "16.5",
    "flavonoids": "16.8",
    "aa_synthesis": "13.1",
    "aa_degradation": "13.2",
    "glycolysis_cystolic": "4.1",
    "glycolysis_plastid": "4.2",
    "glycolysis_other": "4.3"
}

print(PATHWAYS)

def bincode(wildcards):
	bincode = PATHWAYS[wildcards.pathway]
	return bincode

# list trait names
TRAIT = ['ala', 'arg', 'asp', 'gln' , 'glu', 'gly', 'his', 'ile', 'leu', 'lys',
'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'Total', 'ShikFam',
'AspFam', 'AspFam_Asp', 'SerFam', 'GluFam', 'GluFam_glu', 'BCAA', 'PyrFam',
'ala_t', 'arg_t', 'asp_t', 'gln_t', 'glu_t', 'gly_t', 'his_t', 'ile_t', 'leu_t',
'lys_t', 'met_t', 'phe_t', 'pro_t', 'ser_t', 'trp_t', 'tyr_t', 'val_t',
'ile_BCAA', 'leu_BCAA', 'val_BCAA', 'gln_GluFamCorr', 'glu_GluFamCorr',
'his_GluFam', 'pro_GluFamCorr', 'phe_shik', 'trp_shik', 'tyr_shik', 'gly_SerFam',
'ser_SerFam', 'ala_pyr', 'leu_pyr', 'val_pyr', 'asp_AspFam_aspCorr',
'ile_AspFam_aspCorr', 'lys_AspFam_aspCorr', 'met_AspFam_aspCorr', 'thr_AspFam_aspCorr']

# list range of random subsets to include - we used 5000
# RANDOM = list(range(1,5001))
RANDOM = list(range(1,10))

# how many sets of cross validation - we used 5
CV = list(range(1,6))

# fold for each cross validation - we used 10 fold (total of 50 cross validations)
INDEX = list(range(1,11))

# index for null distribution gene group sampling
NULL = np.arange(1, 5001, 50).tolist()

# directory for ldak software
ldak = "../software/ldak5.linux"

################################################################################
# this is a pseudo-rule that collects the target files

rule all:
    input:
        # snp_weightings
        "data/processed/sections/weights.short",
        # pca
        "data/processed/pca.vect",
        # get_cv
        expand("data/processed/cross_validation/cv{cv}.test{index}", \
        cv = CV, index = INDEX),
        expand("data/processed/cross_validation/cv{cv}.train{index}", \
        cv = CV, index = INDEX),
        # gblup_kins
        "models/gblup/kinships.grm.id",
        # gblup_h2
        expand("models/gblup_h2/{trait}.gblup.reml", trait = TRAIT),
        # gblup_reml
        expand("models/gblup/{trait}.cv{cv}.{index}.reml", \
        trait = TRAIT, cv = CV, index = INDEX),
        # gblup_blup
        expand("models/gblup/{trait}.cv{cv}.{index}.profile", \
        cv = CV, index = INDEX, trait = TRAIT),
        # gblup_results
        "reports/gblup.RData",
        # multiblup_pathways
        expand("data/processed/pathways/{pathway}/list1", pathway = PATHWAYS.keys()),
        # multiblup_kins
        expand("data/processed/pathways/{pathway}/partition.list", pathway = PATHWAYS.keys()),
        # multiblup_h2
        expand("models/multiblup_h2/{pathway}/{trait}.multiblup.reml", \
        pathway = PATHWAYS.keys(), trait = TRAIT),
        # multiblup_reml
        expand("models/multiblup/{pathway}/{trait}.cv5.10.reml", \
        pathway = PATHWAYS.keys(), trait = TRAIT),
        # multiblup_blup
        expand("models/multiblup/{pathway}/{trait}.cv5.10.profile", \
        pathway = PATHWAYS, trait = TRAIT),
        # multiblup_results
        "reports/multiblup.RData",
        # null_sampling
        "data/interim/null_group_sizes.txt",
        expand("data/processed/random_sets/null_{null}.txt", null = NULL),
        # calc_kins_control
        expand("data/processed/random_sets/c_{random}/partition.list", random = RANDOM),
        # reml_h2_control
        expand("models/null_h2/c_{random}/{trait}.h2.reml", \
        random = RANDOM, trait = TRAIT)

include: "rules/common.smk"
include: "rules/prep_data.smk"
include: "rules/cross_validation.smk"
include: "rules/gblup.smk"
include: "rules/multiblup.smk"
include: "rules/null_distribution.smk"
