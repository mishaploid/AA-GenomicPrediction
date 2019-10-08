Genomic prediction for free amino acid traits in Arabidopsis seeds
==============================

Data and scripts necessary to run genomic partitioning and prediciton models on free amino acid traits measured in a diverse panel of 313 Arabidopsis lines. 

Software requirements
------------
* R v3.5.1
* PLINK v1.90b4 64-bit
* Miniconda3 (includes conda v4.6.7)
* Snakemake v5.4.2 (install to virtual environment)

To install snakemake in a virtual environment, run:  

`conda env create --name multiblup --file environment.yaml`  

For future use, activate this environment with:
`source activate multiblup`

Data
------------
Edit the config.yaml file to specify the paths for the genotype and phenotype data.

### Genotype data
Before running Snakemake, download the **Arabidopsis Regional Mapping (RegMap) data** ([Horton et al. 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3267885/)):  
`cd data/external`  
`wget https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/250k_snp_data/call_method_75.tar.gz`  
`tar -xvf call_method_75.tar.gz` 

### Phenotype data

#### data/raw/pheno_file
Environment corrected BLUPs for 65 free amino acid traits measured in 313 accessions of _Arabidopsis thaliana_ as reported by [Angelovici et al. 2013](http://www.plantcell.org/content/25/12/4827#sec-12)  

### data/raw/pheno_file_pcs
Phenotypes adjusted for population structure

Snakemake
------------
The snakemake workflow includes a `Snakefile` (specifies steps of the workflow and variables) and rule files in `rules/` to run specific commands. 

### Snakefile 
This file specifies variable names for use in a snakemake workflow. There are notes on how to include specific pathways and MapMan bincodes. Other variables include trait names, the number of random SNP sets, and the number of cross validations to perform.  

The `rule all:` section is a pseudorule that tracks the expected output of each command in the rule files. 

### Rules
rules/common.smk - specifies locaiton of config.yaml file

#### rules/prep_data.smk
- filter and convert genotype data to PED format
- exports TAIR 10 ensembl gene ids for SNP data 
- calculate SNP weightings (these aren't actually used, could skip)
- runs PCA and exports PC adjusted phenotype file

#### rules/cross_validation.smk 
- creates training and testing sets for cross validation

#### rules/gblup.smk
- exports kinship matrix for all SNPs 
- estimates variances and heritability
- genomic prediction and cross-validation (REML + calculating BLUPs)
- summarizes GBLUP output (`reports/gblup.Rdata`)

#### rules/multiblup.smk
- filters exports list of pathway SNPs (includes a 2.5 kb buffer before and after each pathway gene)
- calculates kinship matrices for each pathway (list1) and all remaining SNPs (list2) 
- estimates variances and heritability for each SNP partition
- genomic prediction and cross-validation with multiple kernels/random effects (REML + calculating BLUPs)
- summarizes MultiBLUP output (`reports/multiblup.RData`)

#### rules/null_distribution.smk
- creates 5000 random gene groups with a uniform distribution of SNPs 
- calculates kinship matrices for each random group 
- estimates variances and heritability for each random SNP set 
- summarizes results across all 5000 gene groups (`reports/lr_null_results.csv`)

Notebooks
------------
R notebooks to summarize results. 

#### 01-gblup_results.Rmd
Checks quality of model output and summarizes prediction results for the GBLUP and MultiBLUP models (e.g. proportion of heritability explained, likelihood ratio, prediction accuracy, reliability, bias, MSE)

#### 02-process_null.Rmd
Examines properties of the random gene groups (e.g. distribution of likelihood ratio) and performs quantile regression to establish 95 percentiles for the proportion of heritability explained and likelihood ratio (see nice discussion of this approach in [Edwards et al. 2016](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-015-0132-6))

#### 03-multiblup_results.Rmd
Identifies pathways that pass significance criteria based on comparison to random gene groups with the same number of SNPs (proportion of h<sup>2</sup>, likelihood ratio, and at least a 1% increase in prediction accuracy). 

#### 04-figures.Rmd
Code to create figures used in the manuscript. 


Project Organization (based on Cookiecutter data science)
------------
<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- manuscript documents
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- R notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── environment.yaml   <- The requirements file for reproducing the analysis conda environment
    │
    ├── src                <- Source code for use in this project.
    │   ├── data           <- Scripts to download or generate data
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │                     predictions
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │
    |── Snakefile          <- Snakemake workflow to execute analyses
    │
    └── submit.json        <- Configuration settings to run snakemake on a computing cluster


--------
