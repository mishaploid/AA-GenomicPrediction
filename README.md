Genomic prediction for free amino acid traits in Arabidopsis seeds
==============================

Data and scripts necessary to run genomic partitioning and prediciton models on free amino acid traits measured in a diverse panel of 313 Arabidopsis lines.

Software requirements
------------
* R v3.6.0
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

### Phenotype and covariate data

- **data/raw/aa360_raw_ratios.csv**  
Raw measurements (nmol/mg seed) of 65 free amino acid traits measured in 313 accessions of _Arabidopsis thaliana_ as reported by [Angelovici et al. 2013](http://www.plantcell.org/content/25/12/4827#sec-12)  

- **data/processed/aa360_BLUEs.txt**
Environment adjusted best linear unbiased estimates (BLUEs) for 65 free amino acid traits. Calculated using the `HAPPI-GWAS` pipeline from [Slaten et al. 2020](https://doi.org/10.1093/bioinformatics/btaa589). Check out `notebooks/01-calculate_BLUEs.Rmd` for details.

- **data/processed/aa360_covars.txt**  
Principal components from genotype data to model population structure.

Snakemake
------------
The snakemake workflow includes a `Snakefile` (specifies steps of the workflow and variables) and rule files in `rules/` to run specific commands.

### Snakefile
This file specifies variable names for use in a snakemake workflow. There are notes on how to include specific pathways and MapMan bincodes. Other variables include trait names, the number of random SNP sets, and the number of cross validations to perform.  

The `rule all:` section is a pseudorule that tracks the expected output of each command in the rule files.

To run the workflow, edit cluster configuration settings in `submit.json` then run `submit.sh`

### Rules
rules/common.smk - specifies location of config.yaml file

#### rules/prep_data.smk
- filter and convert genotype data to PED format
- exports TAIR 10 ensembl gene ids for SNP data
- calculate SNP weightings (these aren't actually used, could skip)
- run PCA and exports covariate file (including two PCs here - recommend adjusting depending on data)

#### rules/cross_validation.smk
- create training and testing sets for cross validation
- Note: may need to run this separately on command line, repeating the command via loops/snakemake sometimes does not generate unique folds

#### rules/gblup.smk
- export kinship matrix for all SNPs
- estimate variances and heritability
- genomic prediction and cross-validation (REML + calculating BLUPs)
- summarize GBLUP output (`reports/gblup.Rdata`)

#### rules/multiblup.smk
- filter and export list of pathway SNPs (includes a 2.5 kb buffer before and after each pathway gene)
- calculate kinship matrices for each pathway (list1) and all remaining SNPs (list2)
- estimate variances and heritability for each SNP partition
- genomic prediction and cross-validation with multiple kernels/random effects (REML + calculating BLUPs)
- summarize MultiBLUP output (`reports/multiblup.RData`)

#### rules/null_distribution.smk
- first option: generate 5000 random gene groups with a uniform distribution of SNPs. This is useful if you want to examine influences of partition size on heritability explained/model fit or if you are looking to compare a lot of different partitions with varying size to an empirical distribution (output `reports/lr_null_results.csv`)
- second option: generate 1000 random gene groups for each pathway (excludes pathway SNPs and samples a similar number of SNPs/genes)
- calculate kinship matrices for each random group
- estimate variances and heritability for each random SNP set
- summarize results across all 1000 gene groups (`reports/null_dist_results_pathways/{pathway}_null.csv`)

Notebooks
------------
R notebooks to summarize results.

#### 01-calculate_BLUEs.Rmd


#### 02-gblup_results.Rmd
Checks quality of model output and summarizes prediction results for the GBLUP and MultiBLUP models (e.g. proportion of heritability explained, likelihood ratio, prediction accuracy, reliability, bias, MSE)

#### 03-process_null.Rmd
Examines properties of the random gene groups (e.g. distribution of likelihood ratio) and performs quantile regression to establish 95 percentiles for the proportion of heritability explained and likelihood ratio (see nice discussion of this approach in [Edwards et al. 2016](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-015-0132-6))

#### 04-multiblup_results_by_pathway.Rmd
Identifies pathways that pass significance criteria based on comparison to random gene groups with the same number of SNPs (proportion of h<sup>2</sup>, likelihood ratio, and  increase in prediction accuracy).

#### 05-figures.Rmd
Code to create figures used in the manuscript.


Project Organization (based on Cookiecutter data science)
------------
<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
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
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── environment.yaml   <- The requirements file for reproducing the analysis conda environment
    │
    ├── src                <- Source code for use in this project.
    │   ├── data           <- Scripts to download or generate data
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │                     predictions
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │
    |── Snakefile          <- Snakemake workflow to execute analyses
    │
    └── submit.json        <- Configuration settings to run snakemake on a computing cluster


--------
