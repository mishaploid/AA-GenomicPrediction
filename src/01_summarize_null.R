### AA-GP: process and export summary of null distribution
### Sarah Turner-Hissong
### 07 March 2019

library(tidyverse)
library(quantreg)
library(data.table)

## read in amino acid names + numeric codes 
aa_names <- data.frame(trait = colnames(read.table("data/processed/pheno_file", header = TRUE)[3:56]), id = as.character(1:54))

## import log likelihood for GBLUP model
gblup_llik <- tibble(file = list.files(path = "models/gblup_h2",
                                     pattern = "gblup_\\d+.reml$",
                                     full.names = TRUE)) %>%
  separate(file, sep = "/", into = c("source", "method", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(id = str_replace(id, "gblup_", ""),
         id = str_replace(id, ".reml", ""),
         gblup_llik = as.numeric(X1)) %>%
  left_join(., aa_names, by = "id") %>%
  select(id, trait, gblup_llik)

head(gblup_llik)

## get size of each null gene group
files <- list.files(path = "data/processed/random_sets", pattern = "*.txt", full.names = TRUE)

feature_sizes <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("dir", "dir2", "dir3", "pathway"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  mutate(pathway = str_replace(pathway, "null", "c"),
         pathway = str_replace(pathway, ".txt", "")) %>%
  group_by(pathway) %>%
  summarise(size = length(unique(SNP_ID)),
            n_genes = length(unique(ensembl_gene_id)))
  select(pathway, n_genes, size) %>%
  ungroup() %>%
  write_csv(., "data/processed/random_subsets_summary.csv")

head(feature_sizes)

## list files for reml output for each null gene group
files <- list.files(path = "models/reml_null", pattern = "h2_\\d+.reml$",
                    full.names = TRUE, recursive = TRUE)

## calculate likelihood ratio for null distribution
lr_null <- tibble(file = files) %>%
  separate(file, sep = "/",
           into = c("dir", "source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(id = str_replace(id, "reml_h2_", ""),
         id = str_replace(id, ".reml", ""),
         null_llik = as.numeric(X2)) %>%
  left_join(., aa_names, by = "id") %>%
  left_join(., gblup_llik, by = c("id", "trait")) %>%
  select(trait, id, aa_cat, pathway, Num_Kinships, null_llik, gblup_llik) %>%
  left_join(., null_summary, by = "pathway") %>%
  mutate(group_size = cut(size, breaks = seq(0, 130000, by = 10000),
                          include.lowest = TRUE, dig.lab = 10)) %>%
  mutate(LR = 2 * (null_llik - gblup_llik))

head(lr_null)

## read in variance component estimates
h2_null <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("dir", "source", "pathway", "id"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE, skip = 13)) %>%
  unnest(data) %>%
  mutate(id = str_replace(id, "reml_h2_", ""),
         id = str_extract(id, "\\d+")) %>%
  left_join(., aa_names, by = "id") %>%
  select(trait, pathway, Component, Heritability, Her_SD) 

head(h2_null)

## merge!
null_dist <- left_join(lr_null, h2_null, by = c("trait", "pathway")) %>%
  write_csv(., "data/processed/lr_null_results.csv")

