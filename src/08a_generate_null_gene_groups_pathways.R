################################################################################
## AA-GP: Select SNPs to generate a null distribution
## Author: Sarah Turner Hissong
## Last modified: 9 June 2020

## This script samples gene groups from the data based on annotations from TAIR 10
## All SNPs within a gene are sampled and another gene is sampled until a threshold is reached
## This will be the input for MultiBLUP to generate an empirical null distribution
## Generates 1000 random gene groups with a size comparable to each pathway tested

library(tidyverse)

# use command line arguments
args = commandArgs(trailingOnly = TRUE)

# name of pathway
pathway <- args[1]

# pathway SNPs
snp_ids <- read.table(paste0("data/processed/pathways/", pathway, "/list1"))

# number of SNPs to sample
nsnps <- nrow(snp_ids)

# number of random subsets to export
end <- 1000

cat("pathway:", pathway, "- sampling", end, "random subsets of", nsnps, "snps each\n")

# read in SNP and gene information for all SNPs
snps <- read.table("data/processed/snp_gene_ids_tair10.txt", header = TRUE) %>%
        arrange(chromosome, position) %>%
        filter(!snp_id %in% snp_ids$V1) # remove pathway SNPs

str(snps)

# # read in random uniform distribution for number of SNPs
# snp_number <- read.table("data/interim/gene_groups.txt", header = TRUE)
#
# a function to sample genes randomly up to a specified number of SNPs
# n = numer of SNPs to sample
# bprange = number of basepairs to include before/after a gene

sample_genes <- function(data, nsnps, bprange) {
  out <- data.frame()
  while (nrow(out) < nsnps) {
    # subset snps by genes already selected in output
    set <- na.omit(data[!data$ensembl_gene_id %in% out$ensembl_gene_id,])
    gene <- sample(set$ensembl_gene_id, 1, replace = FALSE) # select a random gene
    snp_sub <- set[set$ensembl_gene_id==gene,] # select all snps within gene
    chr <- snp_sub$chromosome[1]
    mn <- min(snp_sub$position) - bprange
    mx <- max(snp_sub$position) + bprange
    positions <- data[data$position > mn & data$position < mx & data$chromosome==chr,] # select SNPs in gene with X kb buffer
    out <- rbind(out, positions) # combine until threshold number of snps are met
    }
    return(out)
  }

# loop through gene sets specified by commandArgs

for (i in 1:end) {
  if(file.exists(paste0("data/processed/random_sets_pathways/", pathway, "/null_", i, ".txt"))) {
    cat("file already exists for sample", i, "\n")
    next
  }
  tryCatch({

  cat("exporting random gene group", i, "of", end, "\n")

  random <- sample_genes(snps, n = nsnps, bprange = 2500) %>%
  arrange(chromosome, position)

  dir.create(paste0("data/processed/random_sets_pathways/", pathway), showWarnings = FALSE, recursive = TRUE)

  write.table(random, paste0("data/processed/random_sets_pathways/", pathway, "/null_", i, ".txt"))

  # export snplists for LDAK
  dir <- paste0("data/processed/random_sets_pathways/", pathway, "/c_", i)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # SNPs in random gene group
  random %>%
    select(snp_id) %>%
    unique() %>%
    write.table(., paste0(dir, "/list1"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # remaining SNPs
  snps %>%
    select(snp_id) %>%
    filter(!snp_id %in% random$V1) %>%
    write.table(., paste0(dir, "/list2"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  }, error = function(e){cat("ERROR:", conditionMessage(e), "\n")})
}
