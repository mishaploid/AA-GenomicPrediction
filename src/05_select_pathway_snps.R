################################################################################
## AA-GP: Extract gene categories for genomic feature prediction
## Author: Sarah Turner Hissong
## Created: 8 February 2018
## Updated: 18 March 2019

library(tidyverse)
library(fuzzyjoin)

# use command line arguments
args = commandArgs(trailingOnly = TRUE)
print(args[1])
print(args[2])

# read in gene list and annotations (from biomaRt)
gene_list <- read.table("data/processed/gene_list_tair10.txt", header = TRUE)

################################################################################
## Function to extract genes for protein related pathways
## input: data frame with chromosome name, ensembl gene id, external gene name,
## start/end position, and GO category names (name_1006)
################################################################################

extract_genes <- function(genes, category, go_term, snp_data, filename) {
  tmp <- genes[grep(paste0('^', category), genes[,go_term], perl = TRUE),] %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    mutate_at(c("chromosome_name", "ensembl_gene_id"), funs(factor(.))) %>%
    mutate(mn = start_position - 2500,
           mx = end_position + 2500,
           chromosome_name = as.character(chromosome_name))

  snp_data <- snp_data %>% mutate(chromosome = as.character(chromosome))

  snp_tmp <- fuzzy_inner_join(snp_data, tmp, by = c("chromosome" = "chromosome_name",
                                        "position" = "mn",
                                        "position" = "mx"),
                    match_fun = list(`==`, `>=`, `<=`)) %>%
    droplevels() %>%
    dplyr::select(chromosome, snp_id, position, ensembl_gene_id.x,
                  ensembl_gene_id.y, external_gene_name, BINCODE, NAME) %>%
    distinct(., snp_id, .keep_all = TRUE)

  dir.create(paste0("data/processed/pathways/", filename), recursive = TRUE)
  write.table(snp_tmp, paste0("data/processed/pathways/",
                              filename, "/", filename, ".txt"),
              col.names = TRUE, row.names = FALSE)

  return(snp_tmp)
}

# Read in annotations from Mapman

cats <- read.csv("data/external/Ath_AGI_LOCUS_TAIR10_Aug2012.csv", header = TRUE, na.strings = c("", "NA"))
cats$IDENTIFIER <- toupper(cats$IDENTIFIER) # convert locus ids to uppercase
str(cats)
head(cats)

# merge with gene_list
gene_cats <- merge(gene_list, cats, by.x = "ensembl_gene_id", by.y = "IDENTIFIER")

# read in snps and order by Chromosome/Position
snps <- read.table("data/processed/snp_gene_ids_tair10.txt", header = TRUE) %>%
  mutate_at(c("chromosome"), funs(factor(.))) %>%
  mutate(position = as.integer(position)) %>%
  arrange(chromosome, position)

# create directories for pathways
dir <- paste0("data/processed/pathways/", args[1])
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

# pull out annotation categories based on Mapman categories
pathway <- if(args[2] == 0) {
  read.table(paste0(dir, "/", args[1], ".txt"), header = TRUE)
} else {
  extract_genes(gene_cats, c(args[2]), "BINCODE", snps, args[1])
}

head(pathway)

unique(pathway$BINCODE)
unique(pathway$NAME)

# export snplists for LDAK
pathway %>%
  select(snp_id) %>%
  write.table(., paste0(dir, "/list1"), quote = FALSE, row.names = FALSE, col.names = FALSE)

snps %>%
  select(snp_id) %>%
  filter(!snp_id %in% pathway$snp_id) %>%
  write.table(., paste0(dir, "/list2"), quote = FALSE, row.names = FALSE, col.names = FALSE)
