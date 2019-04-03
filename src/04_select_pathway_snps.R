################################################################################
## Extract gene categories for genomic feature prediction
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
  tmp <- genes[grep(paste(category), genes[,go_term]),] %>%
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
                  ensembl_gene_id.y, external_gene_name, BINCODE) %>%
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

# pull out annotation categories based on Mapman categories
pathway <- extract_genes(gene_cats, c(args[2]), "BINCODE", snps, args[1])

# export snplists for LDAK
dir <- paste0("data/processed/pathways/", args[1])
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

pathway %>%
  select(snp_id) %>%
  write.table(., paste0(dir, "/list1"), quote = FALSE, row.names = FALSE, col.names = FALSE)

snps %>%
  select(snp_id) %>%
  filter(!snp_id %in% pathway$snp_id) %>%
  write.table(., paste0(dir, "/list2"), quote = FALSE, row.names = FALSE, col.names = FALSE)


# protein_aa_activation <- extract_genes(gene_cats, c("29.1"), "BINCODE", snps, "protein_activation") # protein.aa activation
# protein_synthesis <- extract_genes(gene_cats, c("29.2"), "BINCODE", snps, "protein_synthesis") # protein.synthesis
# protein_targeting <- extract_genes(gene_cats, c("29.3"), "BINCODE", snps, "protein_targeting") # protein.targeting
# protein_postrans <- extract_genes(gene_cats, c("29.4"), "BINCODE", snps, "protein_posttrans") # protein.postranslational modification
# protein_degradation <- extract_genes(gene_cats, c("29.5"), "BINCODE", snps, "protein_degradation") # protein.degradation
# protein_folding <- extract_genes(gene_cats, c("29.6"), "BINCODE", snps, "protein_folding") # protein.folding
# protein_glyco <- extract_genes(gene_cats, c("29.7"), "BINCODE", snps, "protein_glyco") # protein.glycosylation
# protein_assembly <- extract_genes(gene_cats, c("29.8"), "BINCODE", snps, "protein_assembly") # protein.assembly and cofactor ligation
# aa_transport <- extract_genes(gene_cats, c("34.3"), "BINCODE", snps, "aa_transport")
# tca_cycle <- extract_genes(gene_cats, c("8"), "BINCODE", snps, "tca_cycle")
# degradation_subtilases <- extract_genes(gene_cats, c("29.5.1"), "BINCODE", snps, "degradation_subtilases")
# degradation_autophagy <- extract_genes(gene_cats, c("29.5.2"), "BINCODE", snps, "degradation_autophagy")
# degradation_protease <- extract_genes(gene_cats, c("29.5.3", "29.5.4", "29.5.5", "29.5.7"), "BINCODE", snps, "degradation_proteases")
# degradation_AAA <- extract_genes(gene_cats, c("29.5.9"), "BINCODE", snps, "degradation_AAA")
# degradation_ubiquitin <- extract_genes(gene_cats, c("29.5.11"), "BINCODE", snps, "degradation_ubiquitin")
# e_alt_oxidase <- extract_genes(gene_cats, c("9.4"), "BINCODE", snps, "e_alt_oxidase")
# isoprenoids <- extract_genes(gene_cats, c("16.1"), "BINCODE", snps, "isoprenoids")
# phenylpropanoids <- extract_genes(gene_cats, c("16.2"), "BINCODE", snps, "phenylpropanoids")
# N_containing <- extract_genes(gene_cats, c("16.4"), "BINCODE", snps, "N_containing")
# S_containing <- extract_genes(gene_cats, c("16.5"), "BINCODE", snps, "S_containing")
# flavonoids <- extract_genes(gene_cats, c("16.8"), "BINCODE", snps, "flavonoids")
# aa_synthesis <- extract_genes(gene_cats, c("13.1"), "BINCODE", snps, "aa_synthesis")
# aa_degradation <- extract_genes(gene_cats, c("13.2"), "BINCODE", snps, "aa_degradation")
# # glycolysis <- extract_genes(gene_cats, c("4"), "BINCODE", snps, "glycolysis")
# glycolysis_cystolic <- extract_genes(gene_cats, c("4.1"), "BINCODE", snps, "glycolysis_cystolic") # glycolysis - cystolic branch
# glycolysis_plastid <- extract_genes(gene_cats, c("4.2"), "BINCODE", snps, "glycolysis_plastid") # glycolysis - plastid branch
# glycolysis_other <- extract_genes(gene_cats, c("4.3"), "BINCODE", snps, "glycolysis_other") # glycolysis - unclear or dually targeted
