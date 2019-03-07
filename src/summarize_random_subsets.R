library(tidyverse)

files <- list.files(path = "data/processed/random_sets", pattern = "*.txt", full.names = TRUE)

feature_sizes <- tibble(file = files) %>%
  separate(file, sep = "/", into = c("dir", "dir2", "dir3", "pathway"), remove = FALSE) %>%
  mutate(data = lapply(file, read.table, header = TRUE)) %>%
  unnest(data) %>%
  mutate(pathway = str_replace(pathway, "null", "c"),
         pathway = str_replace(pathway, ".txt", "")) %>%
  group_by(pathway) %>%
  summarise(size = length(unique(SNP_ID)),
            n_genes = length(unique(ensembl_gene_id))) %>%
  select(pathway, n_genes, size) %>%
  ungroup() %>%
  write_csv(., "random_subsets_summary.csv")

