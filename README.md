Genomic prediction for free amino acid traits in Arabidopsis seeds
==============================

Data and scripts necessary to run genomic partitioning and prediciton models on free amino acid traits measured in a diverse panel of 312 Arabidopsis lines. 

Software requirements
------------
* R v3.5.1
* PLINK v1.90b4 64-bit
* Miniconda3 (includes conda v4.6.7)
* Snakemake v5.4.2 (install to virtual environment)

To install snakemake in a virtual environment, run:  

`conda env create --name multiblup --file environment.yaml`  

Setup
------------
Download Arabidopsis Regional Mapping (RegMap) data to `data/external`:  
`cd data/external`  
`wget https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/250k_snp_data/call_method_75.tar.gz`  
`tar -xvf call_method_75.tar.gz`

Use PLINK to remove accessions with missing trait data from genotype data:  

`plink --bfile data/external/plinkGeneOmeSubset --keep data/raw/keep_ids --make-bed --out data/processed/input_nomissing` 


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


