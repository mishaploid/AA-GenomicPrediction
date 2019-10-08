# MultiBLUP workflow for free amino acid traits in Arabidopsis thaliana
# Sarah Turner-Hissong
# Last updated: 18 September 2019

# Phenotype data: Angelovici et al 2013 Plant Cell
# Genotype data: Li et al. 2010; Horton et al. 2012
# Most steps follow MultiBLUP documentation at dougspeed.com/multiblup/

# CAVEATS (in no specific order)
# 1. Full workflow requires a substantial amount of disk space
#    Recommendation is to run just one or a few of the random subsets to start with, then scale up
# 2. Requires phenotype/genotype data to be previously downloaded
# 3. Paths to input files are not generic and will need to be updated if applied to different data

# to execute workflow, use submit.sh file or run the following to include cluster config:
# snakemake --jobs N --rerun-incomplete --latency-wait 60 --cluster-config submit.json
# --cluster "sbatch -p {cluster.p} -o {cluster.o} --cpus-per-task {cluster.cpus-per-task}"

################################################################################
import numpy as np

# define some variables

# dictionary of pathways and MapMan bincodes that will be tested
# to use prespecified files with pathway info, set bincode == 0
# take care when specifying bincodes - want to make sure it is grepping the correct value
# worth double checking when you pull annotations using src/02_select_pathway_snps
# e.g. if you want all annotations in bin 8, just use '8', if you want annotations in '29.5.1',
# specify the end of the line using '29.5.1$' or it will also pull '29.5.11'

# to check this is working correctly in R:
# foo <- tibble(file = files) %>%
# separate(file, sep = "/", into = c("dir", "pathway", "filename"), remove = FALSE) %>%
# filter(!pathway %in% "old") %>% mutate(data = lapply(file, read.table, header = TRUE, colClasses = "character"))
# bar <- foo %>% group_by(pathway) %>% distinct(BINCODE, NAME)

PATHWAYS = {
    "protein_aa_activation": "29.1",
    "protein_synthesis": "29.2",
    "protein_targeting": "29.3",
    "protein_postrans": "29.4",
    "protein_degradation": "(?!29.5.11)29.5",
    "degradation_subtilases": "29.5.1$",
    "degradation_autophagy": "29.5.2$",
    "cys_protease": "29.5.3$",
    "asp_protease": "29.5.4$",
    "ser_protease": "29.5.5$",
    "metallo_protease": "29.5.7$",
    "degradation_AAA": "29.5.9$",
    "degradation_ubiquitin": "29.5.11",
    "protein_folding": "29.6",
    "protein_glyco": "29.7",
    "protein_assembly": "29.8",
    "aa_transport": "34.3$",
    "tca_cycle": "8",
    "e_alt_oxidase": "9.4",
    "isoprenoids": "16.1",
    "phenylpropanoids": "16.2",
    "N_containing": "16.4",
    "S_containing": "16.5",
    "flavonoids": "16.8",
    "aa_synthesis": "13.1",
    "aa_degradation": "13.2",
    "glycolysis": "4",
    "glycolysis_cystolic": "4.1",
    "glycolysis_plastid": "4.2",
    "glycolysis_other": "4.3",
    "aa_syn_nobcat_haplo": "0",
    "aa_deg_bcat_haplo": "0"
}

def bincode(wildcards):
	bincode = PATHWAYS[wildcards.pathway]
	return bincode

# list trait names
TRAIT = ['ala', 'arg', 'asp', 'gln' , 'glu', 'gly', 'his', 'ile', 'leu', 'lys','met',
'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val', 'Total', 'ShikFam', 'AspFam',
'AspFam_Asp', 'SerFam', 'GluFam', 'GluFam_glu', 'BCAA', 'PyrFam', 'ala_t',
'arg_t', 'asp_t', 'gln_t', 'glu_t', 'gly_t', 'his_t', 'ile_t', 'leu_t',
'lys_t', 'met_t', 'phe_t', 'pro_t', 'ser_t', 'trp_t', 'tyr_t', 'val_t',
'ile_BCAA', 'leu_BCAA', 'val_BCAA', 'arg_GluFam', 'gln_GluFam', 'glu_GluFam',
'his_GluFam', 'pro_GluFam', 'phe_ShikFam', 'trp_ShikFam', 'tyr_ShikFam',
'gly_SerFam', 'ser_SerFam', 'ala_PyrFam', 'leu_PyrFam', 'val_PyrFam',
'asp_AspFam', 'ile_AspFam', 'lys_AspFam', 'met_AspFam', 'thr_AspFam']

# list range of random subsets to include - we used 5000
RANDOM = list(range(1,5001))

# index for null distribution gene group sampling
NULL = np.arange(1, 5001, 50).tolist()

# how many sets of cross validation - we used 5
CV = list(range(1,6))

# fold for each cross validation - we used 10 fold (total of 50 cross validations)
INDEX = list(range(1,11))

# path for ldak software
ldak = "../software/ldak5.linux"

################################################################################
# this is a pseudo-rule that collects the target files
# variable names refer to different rules in the workflow
# this is where you can specify wildcard ids for rules

rule all:
    input:
        # linkage disequilibrium adjusted weightings (not used for this workflow)
        snp_weightings = "data/processed/sections/weights.short",
        # principal components analysis to correct for population structure
        pca = "models/pca/ld_pruned.bed",
        pca_result = "data/raw/pheno_file_pcs",
        # create test and training indexes for cross-validation
        get_cv = expand("data/processed/cross_validation/cv{cv}.test{index}", cv = CV, index = INDEX),
        get_cv2 = expand("data/processed/cross_validation/cv{cv}.train{index}", cv = CV, index = INDEX),
        # kinships for GBLUP
        gblup_kins = "models/gblup/kinships.grm.id",
        # heritability estimates for GBLUP
        gblup_h2 = expand("models/gblup_h2/{trait}.gblup.reml", trait = TRAIT),
        # GBLUP cross-validation
        gblup_reml = expand("models/gblup/{trait}.cv{cv}.{index}.reml", trait = TRAIT, cv = CV, index = INDEX),
        gblup_blup = expand("models/gblup/{trait}.cv{cv}.{index}.profile", cv = CV, index = INDEX, trait = TRAIT),
        gblup_results = "reports/gblup.RData",
        # extract pathway information for MultiBLUP
        multiblup_pathways = expand("data/processed/pathways/{pathway}/list1", pathway = PATHWAYS.keys()),
        # kinship matrices for MultiBLUP model
        multiblup_kins = expand("data/processed/pathways/{pathway}/partition.list", pathway = PATHWAYS.keys()),
        # heritability for partitions in MultiBLUP model
        multiblup_h2 = expand("models/multiblup_h2/{pathway}/{trait}.multiblup.reml", pathway = PATHWAYS.keys(), trait = TRAIT),
        # cross-validation for MultiBLUP
        multiblup_reml = expand("models/multiblup/{pathway}/{trait}.cv5.10.reml", pathway = PATHWAYS.keys(), trait = TRAIT),
        multiblup_blup = expand("models/multiblup/{pathway}/{trait}.cv5.10.profile", pathway = PATHWAYS, trait = TRAIT),
        # summary of MultiBLUP results
        multiblup_results = "reports/multiblup.RData",
        # null distribution
        null_sampling = "data/interim/null_group_sizes.txt",
        null_gene_groups = expand("data/processed/random_sets/null_{null}.txt", null = NULL),
        calc_kins_control = expand("data/processed/random_sets/c_{random}/partition.list", random = 5000),
        reml_h2_control = expand("models/null_h2/c_{random}/{trait}.h2.reml", \
        random = 5000, trait = TRAIT),
        null_results = "reports/lr_null_results.csv"

# include rule files with commands to run each step

include: "rules/common.smk"
include: "rules/prep_data.smk"
include: "rules/cross_validation.smk"
include: "rules/gblup.smk"
include: "rules/multiblup.smk"
include: "rules/null_distribution.smk"
