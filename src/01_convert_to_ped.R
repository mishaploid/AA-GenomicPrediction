################################################################################
## AA-GP: convert SNP data to PED/MAP format
## Sarah Turner-Hissong
## Updated 12 March 2019

library(tidyverse)

snp_data <- read.csv("data/external/call_method_75/call_method_75_TAIR9.csv",
                     header = TRUE, skip = 1)

# create the MAP file
# requires Chromosome and Position info
# columns: Chromosome/rsid/physical_distance/Position

map <- snp_data %>%
  select(Chromosome, Positions) %>%
  mutate(rsid = paste0("S", seq_len(n())),
         distance = 0) %>%
  select(Chromosome, rsid, distance, Positions) %>%
  write_delim(., "data/raw/call_method_75_TAIR9.map", col_names = FALSE)

# create the PED file
# rows are individuals
# columns are family/sample/paternal/maternal/sex/phenotype/marker1/marker2/etc

# individual ids
id <- colnames(snp_data[,-c(1:2)])

# get accession ids with phenotype data
acc360 <- read_delim("data/processed/pheno_file", delim = "\t") 

ped <- snp_data[,-c(1:2)] %>%
  # recode genotypes to be diploid
  mutate_if(is.factor, funs(fct_recode(., GG = "G",
                                       TT = "T",
                                       AA = "A",
                                       CC = "C"))) %>%
  t() %>%
  as.data.frame(row.names = FALSE) %>%
  # add columns expected by plink (family/sample id/paternal/maternal/sex/phenotype)
  add_column(., family = "AtRegMap",
             sample = gsub("X", "", id),
             paternal = 0,
             maternal = 0,
             sex = -9,
             phenotype = -9,
             .before = 1) %>%
  # only include accessions with phenotype data
  filter(sample %in% acc360$IID) %>%
  write_delim(., "data/raw/call_method_75_TAIR9.ped", col_names = FALSE)
