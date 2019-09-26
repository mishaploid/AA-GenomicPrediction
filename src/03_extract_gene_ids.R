####################################################################################################
## AA-GP: Extract gene ids for TAIR10 annotation
## Sarah Turner-Hissong
## Created: 12 October 2018
## Last modified: 12 March 2019

# this script extracts gene annotations for Arabidopsis (TAIR 10) using the bioMart
# database and matches annotations to SNPs

# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)
library(data.table)
library(tidyverse)
library(dplyr)

# use command line arguments
args = commandArgs(trailingOnly = TRUE)

# load biomart data for Arabidopsis TAIR 10

listMarts(host = "plants.ensembl.org")
listDatasets(useMart(biomart = "plants_mart", host = "plants.ensembl.org"))
ensembl_plants <- useMart(biomart = "plants_mart", host = "plants.ensembl.org")
aa_genes <- useDataset("athaliana_eg_gene", mart = ensembl_plants)

# load SNP ID and physical position info
snp_query <- read.table(args[1], header = FALSE,
                        col.names = c("chromosome", "snp_id", "physical", "position"))
head(snp_query)
length(snp_query$snp_id)

################################################################################
# extract gene list given SNP locations
gene_list <- unique(getBM(attributes = c("chromosome_name", "ensembl_gene_id", "external_gene_name",
                                         "start_position", "end_position", "name_1006"),
                          filters = c("chromosome_name", "start", "end"),
                          values = list(1:5, snp_query$position, snp_query$position),
                          mart = aa_genes))
head(gene_list)
dim(gene_list)

write.table(gene_list, "data/processed/gene_list_tair10.txt", row.names = FALSE)

################################################################################
## identify SNPs that fall within genic regions
## source: https://stackoverflow.com/questions/24766104/checking-if-value-in-vector-is-in-range-of-values-in-different-length-vector

getValue <- function(x, chr, data) {
  tmp <- data %>%
    filter(start_position <= x, x <= end_position & chromosome_name == chr) %>%
    filter(row_number() == 1)
  return(cbind(position = x, ensembl_gene_id = tmp$ensembl_gene_id))
}

s <- Sys.time()
chr1 <- sapply(X = snp_query$position[snp_query$chromosome==1], FUN = getValue, chr = 1, data = gene_list)
chr2 <- sapply(X = snp_query$position[snp_query$chromosome==2], FUN = getValue, chr = 2, data = gene_list)
chr3 <- sapply(X = snp_query$position[snp_query$chromosome==3], FUN = getValue, chr = 3, data = gene_list)
chr4 <- sapply(X = snp_query$position[snp_query$chromosome==4], FUN = getValue, chr = 4, data = gene_list)
chr5 <- sapply(X = snp_query$position[snp_query$chromosome==5], FUN = getValue, chr = 5, data = gene_list)
e <- Sys.time()
e - s

chr1 <- rbindlist(lapply(chr1, as.data.frame), fill = TRUE)
chr2 <- rbindlist(lapply(chr2, as.data.frame), fill = TRUE)
chr3 <- rbindlist(lapply(chr3, as.data.frame), fill = TRUE)
chr4 <- rbindlist(lapply(chr4, as.data.frame), fill = TRUE)
chr5 <- rbindlist(lapply(chr5, as.data.frame), fill = TRUE)

all <- list(chr1, chr2, chr3, chr4, chr5)

snp_gene_ids <- rbindlist(all, idcol = TRUE)

colnames(snp_gene_ids) <- c("chromosome", "position", "ensembl_gene_id")

snp_gene_ids <- snp_gene_ids %>%
  mutate_at(c("chromosome", "position"), funs(factor(.)))

snp_query <- snp_query %>%
  mutate_at(c("chromosome", "position"), funs(factor(.)))

snps <- merge(snp_gene_ids, snp_query, by = c("chromosome", "position"))

str(snps)

snps$position <- as.numeric(levels(snps$position))[snps$position]

write.table(snps, "data/processed/snp_gene_ids_tair10.txt", row.names = FALSE)
