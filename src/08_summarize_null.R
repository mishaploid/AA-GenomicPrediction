### AA-GP: process and export summary of null distribution
### Sarah Turner-Hissong
### 07 March 2019

library(tidyverse)
library(quantreg)
library(data.table)

# summarize random gene groups (number of SNPs & number of genes)

files <- list.files(path = "data/processed/random_sets", pattern = "*.txt", full.names = TRUE)

feature_sizes <- tibble(file = files) %>%
  separate(file, sep = "/|[.]", into = c("dir", "dir2", "dir3", "pathway", "ext"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  mutate(pathway = str_replace(pathway, "null", "c")) %>%
  group_by(pathway) %>%
  summarise(size = length(unique(snp_id)),
            n_genes = length(unique(ensembl_gene_id))) %>%
  select(pathway, n_genes, size) %>%
  ungroup() %>%
  write_csv(., "reports/random_subsets_summary.csv")

head(feature_sizes)


# import log likelihood for GBLUP model
load("reports/gblup.RData")

## list files for reml output of each random gene group
files <- list.files(path = "models/null_h2", pattern = ".reml$",
                    full.names = TRUE, recursive = TRUE)

## calculate likelihood ratio for null distribution
lr_null <- tibble(file = files) %>%
  separate(file, sep = "/|[.]",
           into = c("dir", "source", "pathway", "trait", "metric", "model"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  filter(Num_Kinships %in% "Alt_Likelihood") %>%
  mutate(null_llik = as.numeric(X2)) %>%
  # correct trait names
  mutate(trait = gsub("Corr", "", trait),
         trait = gsub("_asp", "", trait)) %>%
  left_join(., gblup_llik, by = "trait") %>%
  select(trait, pathway, Num_Kinships, null_llik, gblup_llik) %>%
  left_join(., feature_sizes, by = "pathway") %>%
  mutate(group_size = cut(size, breaks = seq(0, 50000, by = 10000),
                          include.lowest = TRUE, dig.lab = 10),
         LR = 2 * (null_llik - gblup_llik))

head(lr_null)

## read in variance component estimates
files <- list.files(path = "models/null_h2", pattern = ".share$",
                    full.names = TRUE, recursive = TRUE)

h2_null <- tibble(file = files) %>%
  separate(file, sep = "/|[.]", into = c("dir", "source", "pathway", "trait", "metric", "model"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  # correct trait names
  mutate(trait = gsub("Corr", "", trait),
         trait = gsub("_asp", "", trait)) %>%
  left_join(., feature_sizes, by = "pathway") %>%
  select(trait, pathway, size, n_genes, Component, Share, SD)

head(h2_null)

## merge!
null_dist <- left_join(lr_null, h2_null, by = c("trait", "pathway")) %>%
  write_csv(., "reports/lr_null_results.csv")
