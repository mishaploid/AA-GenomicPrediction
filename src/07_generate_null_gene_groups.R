################################################################################
## Select SNPs to generate a null distribution
## Author: Sarah Turner Hissong
## Last modified: 25 December 2018

## This script samples gene groups from the data based on annotations from TAIR 10
## All SNPs within a gene are sampled and another gene is sampled until a threshold is reached
## This will be the input for MultiBLUP to generate an empirical null distribution from 1-50000 markers

library(tidyverse)

# use command line arguments
args = commandArgs(trailingOnly = TRUE)
# end <- as.numeric(args[1]) + 49
end <- 5000
cat("random subsets from ", args[1], "to", end, "\n")

# read in SNP and gene information
snps <- read.table("data/processed/snp_gene_ids_tair10.txt", header = TRUE) %>%
arrange(chromosome, position)

# read in random uniform distribution for number of SNPs
snp_number <- read.table("data/interim/gene_groups.txt", header = TRUE)

# a function to sample genes randomly up to a specified number of SNPs
# n = numer of SNPs to sample
# bprange = number of basepairs to include before/after a gene

sample_genes <- function(data, n, bprange) {
  out <- data.frame()
  while (nrow(out) < n) {
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

for (i in args[1]:end) {
  # if(file.exists(paste0("data/processed/random_sets/null_", i, ".txt"))) {
  #   cat("file already exists for sample", i, "\n")
  #   next
  # }
  # tryCatch({
  # print(snp_number[i,])
  # random <- sample_genes(snps, n = snp_number[i,], bprange = 2500) %>%
  # arrange(chromosome, position)
  # dir.create("data/processed/random_sets", showWarnings = FALSE)
  # write.table(random, paste0("data/processed/random_sets/null_", i, ".txt"))

  # export snplists for LDAK
  dir <- paste0("data/processed/random_sets/c_", i)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  ### edit output
    random <- read.table(paste0(dir, "/list1"))

  # random %>%
  #   select(snp_id) %>%
  #   write.table(., paste0(dir, "/list1"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  snps %>%
    select(snp_id) %>%
    filter(!snp_id %in% random$V1) %>% # filter(!snp_id %in% random$snp_id) %>%
    write.table(., paste0(dir, "/list2"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # }, error = function(e){cat("ERROR:", conditionMessage(e), "\n")})
}
